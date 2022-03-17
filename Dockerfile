FROM continuumio/miniconda3

RUN conda create -n chemodel python==3.7.9

SHELL ["conda", "run", "-n", "chemodel", "/bin/bash", "-c"]


RUN pip3 install --upgrade pip
RUN pip3 install torch==1.8.1 torchvision==0.9.1 torchaudio==0.8.1
RUN pip3 install jupyter
RUN pip3 install pandas==1.2.4
RUN pip3 install joblib==1.1.0
#RUN conda install -c rdkit rdkit==2021.03.1 --override-channels
RUN conda install -c rdkit -c conda-forge rdkit==2021.03.1 --override-channels
RUN pip3 install tqdm==4.60.0
RUN pip3 install scikit-learn==0.24.2
RUN pip3 install pyyaml
RUN pip3 install dgl==0.6.1
RUN pip3 install mordred==1.2.0
RUN pip3 install openpyxl==3.0.5
#RUN pip3 install Django
RUN pip3 install django-pandas==0.6.6
RUN pip3 install django-ckeditor
RUN pip3 install django-import-export
RUN pip3 install django-admin-numeric-filter==0.1.6
RUN pip3 install gunicorn
RUN pip3 install --no-deps django-heroku

RUN pip3 install Django==3.2.12

RUN mkdir /code
WORKDIR /code
ADD . /code

#RUN pip3 list

#CMD python manage.py runserver
CMD gunicorn --bind 0.0.0.0:$PORT config.wsgi
#CMD gunicorn config.wsgi
#CMD python3 manage.py 