from flask import Flask, render_template, request
from werkzeug.utils import secure_filename
from model import Zernike, Defocus
from compute import OpticalPerformance
import os

app = Flask(__name__)

UPLOAD_DIR = 'uploads/'

app.config['UPLOAD_FOLDER'] = UPLOAD_DIR
app.secret_key = 'MySecretKey'

if not os.path.isdir(UPLOAD_DIR):
    os.mkdir(UPLOAD_DIR)

# Allowed file types for file upload
ALLOWED_EXTENSIONS = set(['txt', 'dat', 'npy', 'zer'])

def allowed_file(filename):
    """Does filename have the right extension?"""
    return '.' in filename and \
           filename.rsplit('.', 1)[1] in ALLOWED_EXTENSIONS


@app.route('/upload_demo', methods = ["GET", "POST"])
def index():
    # intitialize the form objects
    form = {
            "intro": Zernike(request.form),
            "defocus": Defocus(request.form)
            }

    result = {
              "intro": None,
              "defocus": None
              }


    filename = None  # default


    ## -----------------------------------------------------------------------------Introduction---------------------------------------------------------------
    if request.method == 'POST' and request.form['btn'] == 'Compute Aberrations':
        # Save uploaded file on server if it exists and is valid
        if request.files:
            file = request.files[form["intro"].filename.name]
            if file and allowed_file(file.filename):
                # Make a valid version of filename for any file ystem
                filename = secure_filename(file.filename)
                file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))

        # call the compute function
        OP = OpticalPerformance(filename=filename)
        result["intro"] = OP.compute_intro(defocus=form["intro"].defocus.data, 
                                  astigmatism=(form["intro"].astigmatism_1.data, form["intro"].astigmatism_2.data),
                                  coma=(form["intro"].coma_1.data, form["intro"].coma_2.data, form["intro"].coma_3.data, form["intro"].coma_4.data),
                                  high_order=form["intro"].higher_order.data)

    ## ----------------------------------------------------------------------------Defocus---------------------------------------------------------------------
    elif request.method == 'POST' and request.form['btn'] == 'Compute Defocused Letters':
        # Save uploaded file on server if it exists and is valid
        if request.files:
            file = request.files[form["defocus"].filename.name]
            if file and allowed_file(file.filename):
                # Make a valid version of filename for any file ystem
                filename = secure_filename(file.filename)
                file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))

        # call the compute function
        OP = OpticalPerformance(filename=filename)
        result["defocus"] = OP.compute_defocus(start = form["defocus"].start.data,
                                    stop = form["defocus"].stop.data,
                                    step = form["defocus"].step.data,
                                    letter_size = form["defocus"].letter_size.data)

    else:
        result = {
          "intro": None,
          "defocus": None
          }


    return render_template("view_upload.html", form=form, result=result)


    


if __name__ == '__main__':
	app.run(debug=True)