# encoding: utf-8
# 2024 R1

# open project
SetScriptVersion(Version="24.1.144")
Open(FilePath="C:/Users/Parrot/Desktop/Recent/Graduate_Projects/Meta_IS_AK_Model/Meta_IS_AK_Test.wbpj")
Extensions.UnloadExtension(
    Id="5e7e01a6-8601-11e8-9f8c-28f10e13ffe6",
    Version="2022.2",
    Format="Binary")
Extensions.UnloadExtension(
    Id="5f463412-bd3e-484b-87e7-cbc0a665e474",
    Version="2024.1",
    Format="Binary")
Extensions.UnloadExtension(
    Id="ba012e44-4f35-4a97-aeff-8fe60efc33c9",
    Version="24.1",
    Format="Binary")
Extensions.UnloadExtension(
    Id="7b0e9e84-396d-4099-9602-2ced9dddc253",
    Version="2024.1",
    Format="Binary")
Extensions.UnloadExtension(
    Id="20180725-3f81-49eb-9f31-41364844c769",
    Version="2024.1",
    Format="Binary")
Extensions.UnloadExtension(
    Id="f3e3da52-fb02-4910-8cc9-980efd047bc6",
    Version="2023.1",
    Format="Binary")

# create new design point
designPoint1 = Parameters.CreateDesignPoint()
parameter1 = Parameters.GetParameter(Name="P1")
designPoint1.SetParameterExpression(
    Parameter=parameter1,
    Expression="DS_FLAP_HEIGHT")    # flap plate thickness 
parameter2 = Parameters.GetParameter(Name="P2")
designPoint1.SetParameterExpression(
    Parameter=parameter2,
    Expression="DS_LO_RIB_WIDTH1")  # LO RIB width1 
parameter3 = Parameters.GetParameter(Name="P3")
designPoint1.SetParameterExpression(
    Parameter=parameter3,
    Expression="DS_LO_RIB_WIDTH2")
parameter4 = Parameters.GetParameter(Name="P4")
designPoint1.SetParameterExpression(
    Parameter=parameter4,
    Expression="DS_LO_RIB_WIDTH3")
parameter5 = Parameters.GetParameter(Name="P5")
designPoint1.SetParameterExpression(
    Parameter=parameter5,
    Expression="DS_LA_RIB_WIDTH")
designPoint1.Retained = True
Parameters.SetBaseDesignPoint(DesignPoint=designPoint1)
backgroundSession1 = UpdateAllDesignPoints(DesignPoints=[designPoint1])
designPoint2 = Parameters.GetDesignPoint(Name="0")
Parameters.SetBaseDesignPoint(DesignPoint=designPoint2)
designPoint1.Retained = False

Parameters.ExportAllDesignPointsData(FilePath="E:/workpack/Matlab/MATLAB_reliability_engineering/Graduate_Projects/Meta_IS_AK_Function/Stress_data.csv")
Save(Overwrite=True)
