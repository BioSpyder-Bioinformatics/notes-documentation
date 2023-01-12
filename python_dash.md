# Dash

Dash is used to create webapp interfaces especially for data visualisation. 
It is composed of *COMPONENTS*, displayed in a page through a *LAYOUT* and interact with each other through *CALLBACK*

- Components
	+ Anything like a menu, button, graph etc
- Layout
	+ Part of dash that allows to customise layout display of components
- Callback
	+ Part that allows components to interact
	+ Composed of a decorator and a function
		* Decorator has output and input
			- Both have the component id, and the property belonging to the component 
	



## Sample page
```py
from dash import Dash, dcc, Output, Input					#pip install dash dash-bootstrap-components
import dash_bootstrap_components as dbc

# Build components
app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTAP])

mytext = dcc.Markdown(children='')
#User input so that whatever the user inputs it gets transfered to the markdown box above. This input has starting text
myinput = dbc.Input(value='# Hello World')

#Customise layout
app.layout = dbc.Container([mytext, myinput])


#Callback component for interactivity
@app.callback(
	Output(mytext, component_property='children'),
	Input(myinput, component_property='value')
	)

def update_title(user_input): #function arguments come from the component property of the Input
	return user_input		  #Returned onkects are assigned to the component property of the Output


#Run app
if __name__ == '__main__':
	app.run_server(port=8051)
```


## Actual sample page, you need to give id and reference by it
```py
from dash import Dash, dcc, Output, Input
import dash_bootstrap_components as dbc

#Components
app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

mytext = dcc.Markdown(children='', id='mytext')
myinput = dbc.Input(value='# This is to have a structure which is working', id='myinput')


#Layout
app.layout = dbc.Container([mytext, myinput])


#Callback
@app.callback(
    Output(component_property='children', component_id='mytext'),
    Input(component_property='value', component_id='myinput')
)


def update_title(user_input):
    return user_input


#Run app
if __name__ == '__main__':
    app.run_server(port=8000)
```




# Add graph component and dropdown menu
```py

from dash import Dash, dcc, Output, Input
import dash_bootstrap_components as dbc

#Graph component
import plotly.express as px

#get data from predefined dataset
df = px.data.medals_long()

# Build components
app = Dash(__name__, external_stylesheets=[dbc.themes.VAPOR])

mytitle = dcc.Markdown(children='Analyse Olympic medals')

mygraph = dcc.Graph(figure={})

dropdown = dcc.Dropdown(options=['Bar Plot', 'Scatter Plot'], value='Bar Plot', clearable=False)#Value is for initial plot displayed


# Customise layout
app.layout = dbc.Container([mytitle, mygraph, dropdown])


# Callback 
@app.callback(
	Output(mygraph, component_property='figure'),
	Input(dropdown, component_property='value')
	)

#Function to swap graphs
def update_graph(user_input):
	if user_input == 'Bar Plot':
		fig = px.bar(data_frame = df, x='nation', y='count', color='medal')

	elif user_input == 'Scatter Plot':
		fig = px.scatter(data_frame=df, x='count', y='nation', color='medal', symbol='medal')

	return fig


#Run app
if __name__ == '__main__':
	app.run_server(port=8051)
```






# Build interactive graph
```py

from dash import Dash, dcc, Output, Input
import dash_bootstrap_components as dbc
import plotly.express as px
import pandas as pd

#Read csv 
df = pd.read_csv('https://github.com')

#Build components
app = Dash(__name__, external_stylesheets=[dbc.themes.LUX])

mytitle = dcc.Markdown(children='')

mygraph = dcc.Graph(figure={})

dropdown = dcc.Dropdown(options=df.columns.values[2:], 
						value='Cesearean Delivery Rate', 
						clearable=False)


# Customise layout
app.layout = dbc.Container([
	dbc.Row([
		dbc.Col([mytitle], width=6)
		], justify='center'),

	dbc.Row([
		dbc.Col([mygraph], width=12)
		]),

	dbc.Row([
		dbc.Col([dropdown], width=6)
		], justify='center')

	], fluid=True)



# Callback
@app.callback(
	Output(mygraph, 'figure'),
	Output(mytitle, 'children'),
	Input(dropdown, 'value')
	)

def update_graph(column_name):
	print(column_name)
	print(type(column_name))

	fig = px.choropleth(data_frame=df,
						locations='STATE',
						locationmode='USA-states',
						scope='usa',
						height=600,
						color=column_name, # give colours based on this value 
						animation_frame='YEAR'
		)
	return fig, '# '+column_name  # We are returning two because we have two outputs, it gets returned in order of declaration




#Run app
if __name__ == '__main__':
	app.run_server(port=8051)
```

























