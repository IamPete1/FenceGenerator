import os
from flask import Flask
from flask import render_template

# Directory of this file
this_path = os.path.dirname(os.path.realpath(__file__))

app = Flask(__name__, static_folder=os.path.join(this_path,'Fences'))

@app.route('/')
def index():
    return render_template('index.html')

if __name__ == "__main__":
    app.run()
