import dash
from dash import Dash, dcc, html, Output, Input, State
from utils.components.sidebar import SideBar

# Initial app declaration, usepages true 
app = Dash(__name__, use_pages=True, pages_folder='utils')

# Set stylesheets
app.config.external_stylesheets = ['./style/base_style.css', './style/extra_style.css']


# App layout
app.layout = html.Div(
    [
        #Div for sidebar (make callback to adapt these 2 styles)
        html.Div(
            SideBar(app),
            id='sidebarDiv'
        ),

        html.Div(
            dash.page_container,
            id='contentDiv'
        )
    ],
    id='mainDiv',
    
    )



# Callbacks






# Run App
if __name__ == '__main__':
    app.run_server(port=9000, debug=True)

























