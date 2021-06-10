'''app

This module contains the flask app.
'''
from flask import Flask
from map.routes import map_bp


app = Flask(__name__)
app.register_blueprint(map_bp)