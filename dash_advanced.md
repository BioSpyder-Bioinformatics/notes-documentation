# Multi-Page App

Folder structure

- app.py
- pages
	+ page1
	+ page2
	+ page3

Page template:
```py
import dash
from dash import dcc, html

dash.register_page(__name__, path='/') #for homepage /, for others


#page layout

layout = html.Div(
		[
			children1,
			children2
		]
	)
```



app.py to allow multi page
```py
import dash
from dash import Dash, dcc, html, Output, Input

# Initial app declaration, usepages true, pages_folder to be specified if named differently than 'pages'; also if 'pages' is in another folder, in this case utils
app = Dash(__name__, use_pages=True, pages_folder='utils')

# Set stylesheets
app.config.external_stylesheets = ['./style/base_style.css', './style/extra_style.css']


# App layout
app.layout = html.Div(
    [
        html.Div('HELLO THERE'),
        #this is for the links, it's a list comprehension to iterate in pre-determined values
        html.Div([
            dcc.Link(page['name']+" | ", href=page['path']) for page in dash.page_registry.values()
        ]),
        #where the app displays the page content 
        dash.page_container
    ]
    )

```





--------------- 
For components:
Pass them the app reference, so they can have their own callbacks!!



--------------
For pages:
Pass tuple (function, address) in description for extra function (soft link them)











