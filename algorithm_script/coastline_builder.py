# -*- coding: utf-8 -*-

from qgis.PyQt.QtCore import QCoreApplication
from qgis.core import (QgsGeometry,
                       QgsProcessing,
                       QgsFeature,
                       QgsLineString,
                       QgsFeatureSink,
                       QgsPoint,
                       QgsProcessingException,
                       QgsProcessingAlgorithm,
                       QgsProcessingParameterFeatureSource,
                       QgsProcessingParameterRasterLayer,
                       QgsProcessingParameterFeatureSink,
                       QgsProcessingParameterNumber,
                       QgsWkbTypes)
from qgis import processing


class CoastlineBuilder(QgsProcessingAlgorithm):
    '''
    Скрипт делает две вещи:
    1. Вычисляет методом линейной интерполяции уровни воды вдоль реки
       на основании известных уровней вверху и внизу по течению.
    2. Генерирует береговые линии на основании линии, построенной 
       на предыдущем шаге и цифровой модели рельефа берегов реки
    '''

    # Задаём имена для исходных слоёв и слоёв результата
    INPUT_LINE = 'INPUT_LINE'
    INPUT_DEM = 'INPUT_DEM'
    UPSTREAM_LEVEL = 'UPSTREAM_LEVEL'
    DOWNSTREAM_LEVEL = 'DOWNSTREAM_LEVEL'
    RIVER_Z_INTERPOLATED = 'RIVER_Z_INTERPOLATED'
    COASTS = 'COASTS'

    def tr(self, string):
        return QCoreApplication.translate('Processing', string)

    def createInstance(self):
        return CoastlineBuilder()

    def name(self):
        return 'coastlinebuilder'

    def displayName(self):
        return self.tr('Построение береговых линий')

    def group(self):
        return self.tr('Гидрография')

    def groupId(self):
        return 'hydrography'

    def shortHelpString(self):
        return self.tr("Построение береговых линий рек с учётом уклона водной поверхности")

    def initAlgorithm(self, config=None):

        self.addParameter(
            QgsProcessingParameterFeatureSource(
                self.INPUT_LINE,
                self.tr('Линия оси водотока'),
                [QgsProcessing.TypeVectorLine]
            )
        )
        
        self.addParameter(
            QgsProcessingParameterRasterLayer(
                self.INPUT_DEM,
                self.tr('Модель рельефа')
            )
        )
        
        upstreamParam = QgsProcessingParameterNumber(
                self.UPSTREAM_LEVEL,
                self.tr('Уровень в начале'),
                QgsProcessingParameterNumber.Double
            )
            
        upstreamParam.setMetadata( {'widget_wrapper':
          { 'decimals': 2 }
        })
        
        self.addParameter(upstreamParam)
        
        downstreamParam = QgsProcessingParameterNumber(
                self.DOWNSTREAM_LEVEL,
                self.tr('Уровень в конце'),
                QgsProcessingParameterNumber.Double
            )
            
        downstreamParam.setMetadata( {'widget_wrapper':
          { 'decimals': 2 }
        })
        
        self.addParameter(downstreamParam)

        self.addParameter(
            QgsProcessingParameterFeatureSink(
                self.RIVER_Z_INTERPOLATED,
                self.tr('Река с интерполированными уровнями')
            )
        )

        self.addParameter(
            QgsProcessingParameterFeatureSink(
                self.COASTS,
                self.tr('Береговые линии')
            )
        )

    def processAlgorithm(self, parameters, context, feedback):

        river = self.parameterAsSource(
            parameters,
            self.INPUT_LINE,
            context
        )
        
        dem = self.parameterAsRasterLayer(
            parameters,
            self.INPUT_DEM,
            context
        )
        
        upstream_level = self.parameterAsDouble(
            parameters,
            self.UPSTREAM_LEVEL,
            context
        )
        
        downstream_level = self.parameterAsDouble(
            parameters,
            self.DOWNSTREAM_LEVEL,
            context
        )

        if river is None:
            raise QgsProcessingException(self.invalidSourceError(parameters, self.INPUT_LINE))
            
        if dem is None:
            raise QgsProcessingException(self.invalidSourceError(parameters, self.INPUT_DEM))

        (river_interpolated, river_interpolated_id) = self.parameterAsSink(
            parameters,
            self.RIVER_Z_INTERPOLATED,
            context,
            river.fields(),
            QgsWkbTypes.LineStringZ,
            river.sourceCrs()
        )

        (coasts, coasts_id) = self.parameterAsSink(
            parameters,
            self.COASTS,
            context,
            river.fields(),
            QgsWkbTypes.LineStringZ,
            river.sourceCrs()
        )
        
        feedback.pushInfo('Проекция {}'.format(river.sourceCrs().authid()))

        if river_interpolated is None:
            raise QgsProcessingException(self.invalidSinkError(parameters, self.RIVER_Z_INTERPOLATED))
        
        if coasts is None:
            raise QgsProcessingException(self.invalidSinkError(parameters, self.COASTS))
                
        # Получаем список объектов исходного слоя
        features = river.getFeatures()
        
        for current, feature in enumerate(features):
            if feedback.isCanceled():
                break
                        
            # Создаём объект оси реки с такой же атрибутивкой, что и у исходника
            z_feature = QgsFeature(feature.fields())
            
            # Создаём объекты левого и правого берегов
            left_coast_feature = QgsFeature()
            right_coast_feature = QgsFeature()
            
            # Извлекаем в переменную длину исходной линии
            length = feature.geometry().length()
            
            # Вычисляем превышение между начальным и конечным уровнями
            common_delta = downstream_level - upstream_level
            
            # Вычисляем коэффициент уклона водной поверхности
            slope = common_delta / length
            
            # Создаём новую линейную геометрию для оси реки
            z_geometry = QgsLineString()
            
            # Создаём геометрию для левого и правого берегов
            left_coast_geometry = QgsLineString()
            right_coast_geometry = QgsLineString()
            
            # Получаем список вершин исходной линии
            points = feature.geometry().asPolyline()
            
            # Вычисляем количество шагов для расчета параметров прогрессбара
            total = 100.0 / len(points) if len(points) > 0 else 0
            
            # Создаём новые точки с координатами такими же, как у исходных вершин,
            # а Z-координату вычисляем линейной интерполяцией.
            # Добавляем созданные точки в новую линейную геометрию.
            total_distance_from_startpoint = 0
            for index, point in enumerate(points):
                
                if index == 0:
                    segment_length = 0
                    azimuth = point.azimuth(points[1])
                elif index == len(points) - 1:
                    segment_length = point.distance(points[index - 1])
                    azimuth = points[-2].azimuth(point)
                else:
                    segment_length = point.distance(points[index - 1])
                    azimuth = (points[index - 1].azimuth(point) + point.azimuth(points[index + 1])) / 2
                
                total_distance_from_startpoint += segment_length
                water_level = (total_distance_from_startpoint * slope) + upstream_level
                z_point = QgsPoint(point.x(), point.y(), water_level)
                z_geometry.addVertex(z_point)
                
                # Получаем размер пикселя на местности по X и по Y
                dem_x_gsd = dem.rasterUnitsPerPixelX()
                dem_y_gsd = dem.rasterUnitsPerPixelY()
                
                # Получаем отметку рельефа в точке на оси реки
                raster_value_at_center, res = dem.dataProvider().sample(point, 1)
                
                # Ищем на растре первый пиксель, отметка которого выше, чем уровень воды в данной точке
                if raster_value_at_center < water_level:
                    left_temp_raster_value = raster_value_at_center
                    step_left = 1
                    while left_temp_raster_value < water_level:
                        left_point_geometry = point.project(dem_x_gsd * step_left, azimuth + 90)
                        left_temp_raster_value, res = dem.dataProvider().sample(left_point_geometry, 1)
                        step_left += 1
                    else:
                        left_point_geometry = z_point.project(dem_x_gsd * step_left, azimuth + 90, 90)
                        left_coast_geometry.addVertex(left_point_geometry)
                    
                    right_temp_raster_value = raster_value_at_center
                    step_right = 1
                    while right_temp_raster_value < water_level:
                        right_point_geometry = point.project(dem_x_gsd * step_right, azimuth - 90)
                        right_temp_raster_value, res = dem.dataProvider().sample(right_point_geometry, 1)
                        step_right += 1
                    else:
                        right_point_geometry = z_point.project(dem_x_gsd * step_right, azimuth - 90, 90)
                        right_coast_geometry.addVertex(right_point_geometry)
                    
                else:
                    feedback.pushInfo(str(False))
                
                # Обновляем прогрессбар
                feedback.setProgress(int(index * total))
            
            # Добавляем в объект оси реки новую линейную геометрию 
            z_feature.setGeometry(z_geometry)
            
            # Добавляем в объекты левого и правого берегов новые геометрии
            left_coast_feature.setGeometry(left_coast_geometry)
            right_coast_feature.setGeometry(right_coast_geometry)
            
            # Добавляем объекты в вывод
            river_interpolated.addFeature(z_feature, QgsFeatureSink.FastInsert)
            coasts.addFeature(left_coast_feature, QgsFeatureSink.FastInsert)
            coasts.addFeature(right_coast_feature, QgsFeatureSink.FastInsert)
        
        return {self.RIVER_Z_INTERPOLATED: river_interpolated, self.COASTS: coasts}
