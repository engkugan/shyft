from ._hbv_stack import *
# Fix up types that we need attached to the model
HbvStateVector.push_back = lambda self, x: self.append(x)
HbvStateVector.size = lambda self: len(self)

HbvModel.cell_t = HbvCellAll
HbvParameter.map_t = HbvParameterMap
HbvModel.parameter_t = HbvParameter
HbvModel.state_t = HbvState
HbvModel.statistics = property(lambda self: HbvCellAllStatistics(self.get_cells()))
HbvModel.hbv_snow_state = property(lambda self: HbvCellHBVSnowStateStatistics(self.get_cells()))
HbvModel.hbv_snow_response = property(lambda self: HbvCellHBVSnowResponseStatistics(self.get_cells()))
HbvModel.priestley_taylor_response = property(lambda self: HbvCellPriestleyTaylorResponseStatistics(self.get_cells()))
HbvModel.hbv_actual_evaptranspiration_response=property(lambda self: HbvCellHBVActualEvapotranspirationResponseStatistics(self.get_cells()))
HbvModel.soil_state = property(lambda self: HbvCellHBVSoilStateStatistics(self.get_cells()))
HbvModel.tank_state = property(lambda self: HbvCellHBVTankStateStatistics(self.get_cells()))
HbvOptModel.cell_t = HbvCellOpt
HbvOptModel.parameter_t = HbvParameter
HbvOptModel.state_t = HbvState
HbvOptModel.statistics = property(lambda self:HbvCellOptStatistics(self.get_cells()))
HbvOptModel.optimizer_t = HbvOptimizer
HbvCellAll.vector_t = HbvCellAllVector
HbvCellOpt.vector_t = HbvCellOptVector
HbvState.vector_t = HbvStateVector
HbvState.serializer_t= HbvStateIo

