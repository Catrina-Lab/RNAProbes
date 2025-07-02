from flask import Flask, render_template, request, redirect, url_for
from sys import argv

app = Flask(__name__)

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/form', methods=['GET', 'POST'])
def form():
    if request.method == 'POST':
        programs = request.form.getlist('programs')
        if not programs:
            # Redirect back to home if no programs were selected. Needed if you navigate directly to it without using index.html first.
            return redirect(url_for('index'))
        return render_template('form.html', programs=programs)
    # If accessed via GET, redirect to home
    return redirect(url_for('index'))

@app.route('/send-request', methods=['POST'])
def submit():
    programs = request.form.get('programs').split(',') #todo: verify values (in the method ig, best place)
    # You can now access the list of selected programs
    # and the other form data to process the request.
    print(programs) # For debugging
    return render_template('request-received.html')

if __name__ == '__main__':
    app.run(debug= (True if len(argv) > 1 else False)) #use debug if any commmand line arguments are inputted
