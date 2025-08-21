#python version >=3.9
from __future__ import annotations

import os

import uuid
from pathlib import Path

from flask import Flask, render_template, request, redirect, url_for, Response, make_response, g, jsonify
from sys import argv

from src.server.DelayedProgram import DelayedProgram
from src.rnaprobes.TFOFinder import tfofinder
from src.rnaprobes.PinMol import pinmol
from src.rnaprobes.smFISH import smFISH
from src.server.program_controller import get_program_object, run_program, set_root
from werkzeug.utils import secure_filename

app = Flask(__name__)
app.config['TEMPLATES_AUTO_RELOAD'] = True

program_names = ["TFOFinder", "PinMol", "smFISH"]
IS_WEB_APP = os.environ.get("IS_WEB_APP")
set_root(Path(__file__).parent)

@app.errorhandler(500)
def application_error(error):
    return str(error), 500

@app.before_request
def track_visitor():
    visitor_id = request.cookies.get("visitor_id")
    g.set_cookie = not visitor_id
    if g.set_cookie:
        visitor_id = str(uuid.uuid4())
    g.visitor_id = visitor_id

@app.after_request
def set_visitor_cookie(response):
    if g.set_cookie:
        response.set_cookie("visitor_id", g.visitor_id, max_age=60*60*24*365*5, httponly = True)  # 5 years
    return response

@app.route('/')
def index():
    return render_template('index.html', is_web_app=IS_WEB_APP, exports=get_exports())

@app.route('/design', methods=['GET', 'POST'])
def design():
    if request.method == 'POST':
        programs = request.form.getlist('programs')
        if not programs:
            # Redirect back to home if no programs were selected. Needed if you navigate directly to it without using index.html first.
            return redirect(url_for('index'))
        return render_template('request.html',
                               programs=programs, all_programs=program_names, is_web_app=IS_WEB_APP, exports=get_exports())
    # If accessed via GET, redirect to home
    return redirect(url_for('index'))

def get_exports():
    return {"TFOFinder": tfofinder.exported_values, "PinMol": pinmol.exported_values, "smFISH": smFISH.exported_values}

@app.route('/send-request', methods=['POST'])
def send_request():
    program = request.args.get("program")
    return run_program(program, g.visitor_id, "The given arguments are invalid",
                         "Something went wrong when calculating your result")

@app.route('/query-result', methods=['GET'])
def query_result():
    program_name = secure_filename(request.args.get("program")) #prevents injection attacks
    id = uuid.UUID(request.args.get("id")) #prevents injection attacks
    program = get_program_object(program_name)
    output_dir = program.output_dir / str(id)
    if isinstance(program, DelayedProgram):
        return program.get_current_result(output_dir, id)
    else:
        return f"Can only query a DelayedProgram. {program_name} is not a DelayedProgram", 400

@app.route('/legal')
def legal():
    return render_template('legal.html')

@app.route('/about')
def about():
    return render_template('about.html')

@app.route('/contact')
def contact():
    return render_template('contact.html')

if __name__ == '__main__':
    app.run(debug= (True if "-d" in argv or "--debug" in argv else False)) #use debug if any commmand line arguments are inputted
