import dash
from dash import Dash, dcc, html, Output, Input
# from utils.components.sidebar import SideBar
# from utils.pages.homepage import HomePage

# Initial app declaration, usepages true 
app = Dash(__name__, use_pages=True, pages_folder='utils')# 

# Set stylesheets
app.config.external_stylesheets = ['./style/base_style.css', './style/extra_style.css']


# App layout
app.layout = html.Div(
    [
        html.Div('HELLO THERE'),
        # SideBar('This is the sidebar'), 
        # HomePage('This is homepage'),
        #this is for the links, it's a list comprehension to iterate in pre-determined values
        html.Div([
            dcc.Link(page['name']+" | ", href=page['path']) for page in dash.page_registry.values()
        ]),
        #where the app displays the page content 
        dash.page_container
    ],
    #id='mainDiv',
    
    )



# Callbacks




# Run App
if __name__ == '__main__':
    app.run_server(port=9000, debug=True)

























