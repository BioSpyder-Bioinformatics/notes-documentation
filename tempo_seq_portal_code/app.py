from dash import Dash, dcc, html, Output, Input



# Initial app declaration
app = Dash(__name__)

# Set stylesheets
app.config.external_stylesheets = ['./style/base_style.css', './style/extra_style.css']

# App layout
app.layout = html.Div('HELLO THERE')



# Callbacks




# Run App
if __name__ == '__main__':
    app.run_server(port=9000, debug=True)

























