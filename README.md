# Ionic conductivity predictor 
- Database of Li conducting liquid and solid electrolytes (organic)
    - polymers, liquids, etc
- Predict thier ionic conductivity by machine learning
- Related paper
    - https://pubs.acs.org/doi/10.1021/jacs.9b11442
    - Database and models are slightly different from the original paper
        - mainly because of copyright and available module version issues

# How to use?
- [DEMO server](https://ionpred.herokuapp.com/admin/) is available! (just wait for ca. 30 seconds to load the page)
    - IMPORTANT NOTES about the server
        - Login info
            - username: user
            - pass: user
        - Database will reset every day
        - Response is slow because it runs on a free server


# Installization
## For computer beginners
- The system works on [Django framework](https://docs.djangoproject.com) of Python
    - Basic knowledge of Python and Django would be needed to run the program

### On your server 
1. Clone this repositry
2. Setup Python environment according to "requirements.txt"
    - Or, manually run the commands [here](misc/conda_command) 
3. Run server
    - ```python manage.py runserver```
    - Or, by other command, such as 
        - ```gunicorn -b :8765 config.wsgi```
4. Access website
5. You can login the site with
    - Username: user
    - Pass: user

### On [Heroku](https://heroku.com/) (via docker)
1. Clone this repositry
2. Login heroku via CLI
3. Run the following commands
    - ```heroku create [your heroku project name]```
    - ```heroku container:login```
    - ```heroku stack:set container```
    - ```heroku container:push web -a [your heroku project name]```
    - ```heroku container:release web -a [your heroku project name]```

### Docker (not checked)
1. Clone this repositry
2. Build image
    - ```docker build -t ion .```
3. Run (e.g., @ PORT=8000)
    - ```docker run -e PORT=8000 ion```

# Version
- 2022.3.17 First prototype

# Author
- Kan Hatakeyama-Sato
- Waseda University
- https://kanhatakeyama.github.io/
- satokan@toki.waseda.(japan)