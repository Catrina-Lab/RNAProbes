from flask import Flask, render_template, request

app = Flask(__name__)

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/form/<program>')
def form(program):
    return render_template('form.html', program=program)

@app.route('/submit', methods=['POST'])
def submit():
    # This is where you'll handle the form submission
    # and call your backend scripts.
    # For now, we'll just display a confirmation page.
    return render_template('request-received.html')

if __name__ == '__main__':
    app.run(debug=True)
