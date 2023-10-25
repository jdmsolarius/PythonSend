from os import environ
from ProteinViewer import app  # Import the 'app' object from your app module

if __name__ == '__main__':
    HOST = environ.get('SERVER_HOST', 'localhost')
    
    try:
        PORT = int(environ.get('SERVER_PORT', '5025'))
    except ValueError:
        PORT = 5025

    app.run(threaded=False, host='0.0.0.0', port=PORT)