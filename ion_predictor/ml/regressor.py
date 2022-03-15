#model definition
from sklearn.ensemble import RandomForestRegressor
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder,StandardScaler
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import make_pipeline


numeric_preprocessor = Pipeline(
    steps=[
        ('imputer', SimpleImputer(strategy='most_frequent')),
        ('scaler', StandardScaler()),
    ]
)

categorical_preprocessor = Pipeline(
    steps=[
        ('imputer', SimpleImputer(strategy='most_frequent')),
        ('onehot', OneHotEncoder(handle_unknown='ignore')),
    ]
)


def initiate_regressor(number_cols,category_cols):
    pipeline = make_pipeline(
        ColumnTransformer(
            [
                ('numerical', numeric_preprocessor, number_cols),
                ('categorical', categorical_preprocessor, category_cols)
            ]
        ),

        RandomForestRegressor( n_jobs = -1,max_depth=20)
    )

    return pipeline