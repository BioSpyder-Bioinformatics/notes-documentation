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

# Login route - remember to add methods
@app.route('/login', methods=['GET', 'POST'])
def login():
	form = LoginForm()

	# Check if the form is valid or not
	if form.validate_on_submit():
		return render_template('login.html', form=form)

# Look below for template of login.html, he uses variables and passes the login/register forms from there

# Register route
@app.route('/register', methods=['GET', 'POST'])
def register():
	form = RegisterForm()

	if form.validate_on_submit():
		return render_template('register.html', form=form)



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

<!-- Make form block -->
<form action="{{ url_for('login') }}" method="POST">
	<!-- Protects against some sort of cyberattack related to form submission -->
	{{ form.hidden_tag() }}
	{{ wtf.form_field(form.email) }}
	{{ wtf.form_field(form.password) }}
	<input type="submit" value="login">
</form>
<!-- Button to submit -->

{% endblock %}

```


(For reference register.html)
```html
{% extends 'bootstrap/base.html' %} <!-- This is for the bootstrap theme -->
{% import 'bootstrap/wtf.html' as wtf %} <!-- For the input fields -->

<!-- Title of the page -->
{% block title %} Register {% endblock %}

<!-- Block for page itself -->
{% block content %}

<!-- Make form block -->
<form action="{{ url_for('register') }}" method="POST">
	<!-- Protects against some sort of cyberattack related to form submission -->
	{{ form.hidden_tag() }}
	{{ wtf.form_field(form.email) }}
	{{ wtf.form_field(form.password) }}
	{{ wtf.form_field(form.repeat_pass) }}
	<input type="submit" value="login">
</form>
<!-- Button to submit -->

{% endblock %}

```





## Manage connection to SQL DB
`pipenv install flask-sqlalchemy flask-migrate`

Again in app.py
```py
from flask import Flask, render_template
from flask_wtf import FlaskForm
from wtforms import StringField, PasswordField
from wtforms.validators import Length, Email

# New
from flask_sqlalchemy import SQLAlchemy
from flask_migrate import Migrate

app = Flask(__name__)
app.config['SECRET_KEY'] = 'this_is_a_secret_key'

#-----

# Give address of db (SECOND STEP)
app.config["SQLALCHEMY_DATABASE_URI"] = "sqlite:///sqlite.db"

# Configure DB on app
db = SQLAlchemy(app)
migrate = Migrate(app, db)

# Make class for our User (inherits from db)
class User(db.Model):
	email = db.Column(db.String(128), primary_key=True)
	password = db.Column(db.String(128))
#-----


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

# Login route - remember to add methods
@app.route('/login', methods=['GET', 'POST'])
def login():
	form = LoginForm()

	# Check if the form is valid or not
	if form.validate_on_submit():
		return render_template('login.html', form=form)

# Look below for template of login.html, he uses variables and passes the login/register forms from there

# Register route
@app.route('/register', methods=['GET', 'POST'])
def register():
	form = RegisterForm()

	if form.validate_on_submit():
		return render_template('register.html', form=form)



if __main__ == '__name__':
	app.run()
 ```


Now to connect our app to DB
- We want to have the pipenv activated `pipenv shell`
- `flask db init`
	+ This adds the 'migrations' folder, which manages all the interactions
- Now the DB is connected
- To create a new connection version (you might want to add things at different stages)
	+  `flask db migrate -m "initial"`
	+ Created a new connection named initial
	+ This will override old versions


Now that the migrations are set, we need to indicate where the DB is. Add in app.py before declaring the db variable:
`app.config["SQLALCHEMY_DATABASE_URI"] = "sqlite:///sqlite.db"` 
(no need to declare the file beforehand)
Then on terminal
`flask db upgrade`





### Add Login with flask-login
`pipenv install flask-login`

Looking at app.py - added changes are marked with a comment!

```py
# User class +++++++ We need to add an id bc is required by flasklogin
class User(db.Model):
	#                  NEW
	id=db.Column(db.Integer, primary_key=True) # Need the autoincrement!
	email = db.Column(db.String(128), nullable=False) # NO MORE PRIMARY KEY. but not null
	password = db.Column(db.String(128), nullable=False)
 ```

*Now you can* `flask db migrate -m "message"` in the terminal to update the tables (only the user one in this case, might want to add the id to start with!) (Deffo from beginning, you'd need to change everything in the automatically done fns)





```py
from flask import Flask, render_template, redirect     #!!!! new
from flask_wtf import FlaskForm
from wtforms import StringField, PasswordField
from wtforms.validators import Length, Email
from flask_sqlalchemy import SQLAlchemy
from flask_migrate import Migrate
# NEW
from flask.helpers import url_for
from flask_login import LoginManager, login_user, logout_user
from flask_login.mixins import UserMixin
from werkzeug.security import generate_password_hash, check_password_hash


app = Flask(__name__)
app.config['SECRET_KEY'] = 'this_is_a_secret_key'
app.config["SQLALCHEMY_DATABASE_URI"] = "sqlite:///sqlite.db"
db = SQLAlchemy(app)
migrate = Migrate(app, db)


####
#set up login manager
login = LoginManager()
login.init_app(app)


#Add user loader fn!! In charge of retrieving users, given the id (Called by flask not us)
@login.user_loader
def user_loader(user_id):
	return User.query.filter_by(id=user_id).first()

# ++ need to add methods other to user class, so extend it with UserMixin (instead of add manually)

#### 


# User class 
class User(db.Model, UserMixin): ###############
	id=db.Column(db.Integer, primary_key=True)
	email = db.Column(db.String(128), nullable=False)
	password = db.Column(db.String(128), nullable=False)


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



# Login route  ################################ NEW
@app.route('/login', methods=['GET', 'POST'])
def login():
	form = LoginForm()

	# Check if the form is valid or not
	if form.validate_on_submit():
		# If form is valid, we want to query the table by email, and check the hashed pass against the obj
		user = User.query.filter_by(email=form.email.data).first()
		
		# If this is true the check passed
		if check_password_hash(user.password, form.password.data):
			login_user(user)

			# Redirect to homepage, in this case the app, check how's best!
			return redirect(url_for('index')) #########################



	return render_template('login.html', form=form)


# Register               #################################################
@app.route('/register', methods=['GET', 'POST']) ######### NEW
def register():
	form = RegisterForm()

	# Check form validation AND repeated passwords are the same
	if form.validate_on_submit() and form.password.data == form.repeat_password.data:


		# Make user model to return (with data of forms)
		user = User(
			email = form.email.data,
			# Save password in hashed form
			password = generate_password_hash(form.password.data)
			)

		# Save our user to db
		db.sesssion.add(user)
		db.session.commit()

		return redirect(url_for('login'))        ################# TO be honest I would redirect to login after registering! (was index) 

	return render_template('register.html', form=form)



# New ENDPOINT
# Endpoint for logout (cleans cookies)
@app.route("/logout")
def logout():
	logout_user()
	return redirect(url_for('login'))



if __main__ == '__name__':
	app.run()
 ```







## Integrate Flask and Dash
(Classic installation eg pip install dash pandas)

- Create a new folder for the dash application!
	+ Make a `__init__.py` file (its going to be the old app.py, which now is the flask app)

```py

# All the imports and stuff
# TO PROTECT ROUTES
from flask_login.utils import login_required


# Create a function that receives a flask application and will connect dash and flask!

def create_dash_application(flask_app):
	dash_app = dash.Da sh(
			# Declare on what server it runs
			server=flask_app,
			name='name',
			# Give base url (so you can protect it)
			url_base_pathname='/dash/' # Change this to maybe app!
		)
	
	dash_app.layout = [...]


	#possibly callbacks?






	# HERE PROTECT ROUTES!!!
	for view_function in dash_app.server.view_functions:
		# If it has 'dash' in front of it
		if view_function.startswith(dash_app.config_url_base_pathname):
			# Set it so to have to be logged in!
			dash_app.server.view_functions[view_function] = login_required(dash_app.server.view_functions[view_function])

	return dash_app




```



IN APP.PY
```py
# All imports....
from dash_application import create_dash_application # dash application is the folder



# After all the initialisation, before user loader
create_dash_application(app)




```







































