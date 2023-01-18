from dash import Dash, dcc, html, Output, Input
from utils.components.sidebar import SideBar
from utils.pages.homepage import HomePage

# Initial app declaration
app = Dash(__name__)

# Set stylesheets
app.config.external_stylesheets = ['./style/base_style.css', './style/extra_style.css']

# App layout
app.layout = html.Div(
    [SideBar('This is the sidebar'), 
    HomePage('This is homepage')],
    id='mainDiv'
    )



# Callbacks




# Run App
if __name__ == '__main__':
    app.run_server(port=9000, debug=True)

























