from flask import (
    Blueprint, flash, redirect, render_template, request, url_for, session, g
)
from werkzeug.security import generate_password_hash, check_password_hash
from flask_app.db import get_db
import functools


bp = Blueprint('index', __name__, url_prefix='/')

@bp.before_app_request
def load_logged_in_user():
    user_id = session.get('user_id')

    if user_id is None:
        g.user = None
    else:
        g.user = get_db().execute(
            'SELECT * FROM user WHERE user_id = ?', (user_id,)
        ).fetchone()
        
def login_required(view):
    @functools.wraps(view) 
    def wrapped_view(**kwargs):
        if g.user is None:
            return redirect(url_for('index.login'))
        
        return view(**kwargs)
    
    return wrapped_view

@bp.route('/')
def index():
    return render_template('index/index.html')

@bp.route('/login', methods = ['GET', 'POST'])
def login():
    if request.method == 'POST':
        username = request.form['username']
        password = request.form['password']
        db = get_db()
        error = None
        user = db.execute(
            'SELECT * FROM user WHERE username = ?', (username,)
        ).fetchone()
        
        if user is None:
            error = 'Incorrect username.'
        elif not check_password_hash(user['password'], password):
            error = 'Incorrect password.'

        if error is None:
            session.clear()
            session['user_id'] = user['user_id']
            return redirect(url_for('index.home'))

        flash(error)

    return render_template('index/login.html')

@bp.route('/logout')
def logout():
    session.clear()
    return redirect(url_for('index.index'))


@bp.route('/home')
@login_required
def home():
    return render_template('index/home.html')



