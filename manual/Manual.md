# TOC
- [TOC](#toc)
- [1. Installization](#1-installization)
  - [NOTE: For programming beginners](#note-for-programming-beginners)
  - [Choice A: On your server](#choice-a-on-your-server)
  - [Choive A': On your server via Docker (not checked)](#choive-a-on-your-server-via-docker-not-checked)
  - [Choice B: On Heroku (remote webserver via docker)](#choice-b-on-heroku-remote-webserver-via-docker)
- [2. Overview for using the program](#2-overview-for-using-the-program)
  - [GUI](#gui)
  - [Jupyter](#jupyter)
- [3. Detailed steps](#3-detailed-steps)
  - [3.1 Launch server](#31-launch-server)
- [4. TODO & issues](#4-todo--issues)



# 1. Installization
## NOTE: For programming beginners
- The system works on [Django framework](https://docs.djangoproject.com) of Python
    - Basic knowledge of Python and Django would be needed to run the program
- However, you can try [VirtualBox](https://www.virtualbox.org/).
    - [Virtual Box image](https://drive.google.com/drive/folders/1blh2ysu-766BYBRP9J_iByIbmjkD4GW-?usp=sharing) (2022/4/5, ca. 10 GB)
    - Files were not checked very carefully
    - It may not be the newest version

## Choice A: On your server 
1. Clone this repositry
    1. For instance,
        - ```gh repo clone KanHatakeyama/ion_predictor```
    2. Unzip database
        - ```7z x db.7z```
2. Setup Python environment according to "requirements.txt"
    - Or, manually run the commands [here](../misc/conda_command) 
3. Run server
    - ```python manage.py runserver```
    - Or, by other command, such as 
        - ```gunicorn -b :8765 config.wsgi```
4. Access website
5. You can login the site with
    - Username: user
    - Pass: user

## Choive A': On your server via Docker (not checked)
1. Clone this repositry
    1. For instance,
        - ```gh repo clone KanHatakeyama/ion_predictor```
    2. Unzip database
        - ```7z x db.7z```
2. Build image
    - ```docker build -t ion .```
3. Run (e.g., @ PORT=8000)
    - ```docker run -e PORT=8000 ion```


## Choice B: On [Heroku](https://heroku.com/) (remote webserver via docker)
1. Clone this repositry
    1. For instance,
        - ```gh repo clone KanHatakeyama/ion_predictor```
    2. Unzip database
        - ```7z x db.7z```
2. Login heroku via CLI
    - ```heroku login --interactive ```
3. Run the following commands
    - ```heroku create [your heroku project name]```
    - ```heroku container:login```
    - ```heroku stack:set container```
    - ```heroku container:push web -a [your heroku project name]```
    - ```heroku container:release web -a [your heroku project name]```

# 2. Overview for using the program 
## GUI
- Run server
- Access the program
    - e.g., http://127.0.0.1:8000/
- Edit chemical and electrolyte data
    - You can import and export data as e.g., xlsx and csv

## Jupyter
- Launch [jupyter notebook](../prepare_model.ipynb)
- You can tune neural descriptors, etc


# 3. Detailed steps
## 3.1 Launch server
- 1. Launch the server
    - Example commands
      - Normal launch
        - ```python manage.py runserver```
      - For remote access with the port of 8765
          - ```python manage.py runserver 0.0.0.0:8765```
    - Please check [Django document](https://docs.djangoproject.com/en/4.0/ref/django-admin/) for details about launching the command
        - Note that this command is only for deveplment

- 2. Access the server with your browser
    - e.g., http://127.0.0.1:8000/
    - You will reach a login page


(under devlopment)

# 4. TODO & issues
- Error occurs during exporting/importing large records (e.g., dump all composite data)
    - this seems to be induced by a timeout error of wsgi
    - launching server by django will not cause the error
        - ```python manage.py runserver```

