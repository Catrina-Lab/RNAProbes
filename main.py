from flask import Flask, render_template, request, redirect, url_for
from sys import argv
from src.TFOFinder import tfofinder


app = Flask(__name__)
program_dict = { #get args, validate args, return value
    'tfofinder': (lambda req: (request.form.get("tfo-probe-length")),
                tfofinder.validate_arguments,
                  tfofinder.calculate_result), #temp functions
    'smfish': (lambda x: x, lambda x: x),
    'pinmol': (lambda x: x, lambda x: x),
}

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
        return render_template('request_form.html', programs=programs)
    # If accessed via GET, redirect to home
    return redirect(url_for('index'))

@app.route('/send-request', methods=['POST'])
def send_request():
    programs = request.form.get('programs').split(',')
    # You can now access the list of selected programs
    # and the other form data to process the request.
    print(request.form) # For debugging
    file = request.files.get("ct-file")

    for prog_name in programs:
        args = program_dict[prog_name][0](request)
        if program_dict[prog_name][1](*args): program_dict[prog_name][2](*args)
        else : pass #todo: send fail message
    return render_template('request-received.html')

if __name__ == '__main__':
    app.run(debug= (True if len(argv) > 1 else False)) #use debug if any commmand line arguments are inputted
