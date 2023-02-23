# Integrate your Dash application with Flask

### Structure of Flask app

*Required packages*
- flask (comes with dash)
- flask-wtf (library for submitting forms)
- email-validator (self explanatory)

In top level (app.py)
```py
from flask import Flask, render_template
# This underneath is required for standard authentication of forms
from flask_wtf import FlaskForm
from wtforms import StringField, PasswordField
from wtforms.validators import Length, Email

app = Flask(__name__)
# secret key for the app
app.config['SECRET_KEY'] = 'this_is_a_secret_key'



# Class for logging in
class LoginForm(FlaskForm):
	email = StringField('email', validators=[Email()])
	password = PasswordField('password', validators=[Length(min=5)])

#  Class for registering a user
class RegisterForm(FlaskForm):
	email = StringField('email', validators=[Email()])
	password = PasswordField('password', validators=[Length(min=5)])
	repeat_password = PasswordField('repeated_password', validators=[Length(min=5)])





# Declare a route for the server, returning a html template to render
@app.route('/')
def index():
	return render_template('index.html')


@app.route('/login')
def login():
	return
# Look below for template of login.html, he uses variables and passes the login/register forms from there






if __main__ == '__name__':
	app.run()

```
Run app:
`flask run`

(index.html is any html file)



(For reference login.html)
```html
{% extends 'bootstrap/base.html' %} <!-- This is for the bootstrap theme -->
{% import 'bootstrap/wtf.html' as wtf %} <!-- For the input fields -->

<!-- Title of the page -->
{% block title %} Login {% endblock %}

<!-- Block for page itself -->
{% block content %}

<!-- Protects against some sort of cyberattack related to form submission -->
{{ form.hidden_tag() }}

{{ wtf.form_field(form.email) }}
{{}}



{% endblock %}

```





















































