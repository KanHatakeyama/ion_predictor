#from .ion_predictor.ml import pretrain_descriptors 
#from .ion_predictor.ml.auto_trainer import auto_prepare_model
#from .ion_predictor.ml.NeuralDescriptor import NeuralDescriptor
from ..composite.auto_data_preparer import load_ion_excel
import yaml
import joblib
import numpy as np
import copy

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
