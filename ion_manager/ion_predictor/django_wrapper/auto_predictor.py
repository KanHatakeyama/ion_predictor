#from .ion_predictor.ml import pretrain_descriptors 
#from .ion_predictor.ml.auto_trainer import auto_prepare_model
#from .ion_predictor.ml.NeuralDescriptor import NeuralDescriptor
from ..composite.auto_data_preparer import load_ion_excel
import yaml
import joblib
import numpy as np
import copy

def pre_convert(composite_df,compound_df):
    # convert original django-type df into standard format

    # prepare "composition" column (e.g., s1/s2/s4)
    component_data=list(composite_df[[f"component{i+1}" for i in range(6)]].values)

    composition_data=[]
    for component_record in component_data:
        component_record=[i for i in component_record if i is not None]
        composition="/".join(component_record)
        composition_data.append(composition)

    composite_df["composition"]=composition_data

    #rename columns
    compound_df=compound_df.rename(columns={'title': 'ID'})


    composite_df=composite_df.rename(columns={
                                "title":"ID",
                                'inorg_contain_ratio': "inorg_contain_ratio(wt)",
                                "temperature":"Temperature",
                                "conductivity":"Conductivity"
                                })

    compound_df=compound_df.replace([None], [np.nan])
    composite_df=composite_df.replace([None], [np.nan])

    return composite_df,compound_df



#predict conductivity from dataframe parsed from django admin
def predict(composite_df,compound_df):

    composite_df,compound_df=pre_convert(composite_df,compound_df)

    setting_path="settings.yaml"
    with open(setting_path) as file:
        settings= yaml.safe_load(file)

    #load regressor
    model,X_columns=joblib.load(settings["regressor_path"])

    #prepare machine learnable df
    test_df=load_ion_excel(settings,compound_df=compound_df,composite_df=composite_df)
    test_X=test_df.drop([y_label,"ID"],axis=1)

    #add lacking columns
    lacking_columns=set(X_columns)-set(test_X.columns)
    for c in lacking_columns:
        test_X[c]=np.nan
    test_X=test_X.sort_index(axis=1, ascending=False)


    pred_y=model.predict(test_X)
    test_df["predict"]=pred_y

    return test_df


###################

setting_path="settings.yaml"

with open(setting_path) as file:
    settings= yaml.safe_load(file)


y_label=settings["y_label"]

def screen_predict(path):

    model,X_columns=joblib.load(settings["regressor_path"])

    new_settings=copy.deepcopy(settings)
    new_settings["excel_path"]=path

    test_df=load_ion_excel(new_settings,composite_sheet_name="DatabaseForTest")
    test_X=test_df.drop([y_label,"ID"],axis=1)
    #test_y=test_df[[y_label]]

    lacking_columns=set(X_columns)-set(test_X.columns)
    for c in lacking_columns:
        test_X[c]=np.nan
    test_X=test_X.sort_index(axis=1, ascending=False)


    pred_y=model.predict(test_X)

    test_df["predict"]=pred_y

    return test_df
