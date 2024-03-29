# VExD -- the Virus Expression Database

This project contains the source code for the [VExD website](https://vexd.cchmc.org).  It does not include any data or tools for creating the data.

This site is written using the [Flask](https://palletsprojects.com/p/flask/) framework, which means that it can be run in any [WSGI](https://en.wikipedia.org/wiki/Web_Server_Gateway_Interface) environment.

## Prerequisites

* An accessible Mongo database, with a collection called `vexd`.
* Python 3.8+
* The python packages Flask, marshmallow, webargs, pymongo, matplotlib, and roman (these will be installed by pip).
* Pip will also install the pandas and pyarrow packages, used for the create-downloads command

## Installing

Clone the repository and run the standard Python install (I **always** work in virtual environments):

```
$ git clone https://github.com/pdexheimer/vexd/vexd.git
$ cd vexd
$ python -mvenv env
$ source env/bin/activate
(env) $ pip install -e .
```

## Configuring

VExD will parse a file in `instance/config.py` for runtime configuration.  The format of this file is just a series of `key = value` lines.  The most important variables are:

* `SECRET_KEY` - Used [by Flask](https://flask.palletsprojects.com/en/2.0.x/config/?highlight=secret_key#SECRET_KEY) as a key for signing secure cookies and the like. VExD doesn't use cookies, but it's still good practice to set this properly. If this is the only variable you need, you can generate the file by running `python -c 'import os; print("SECRET_KEY =", os.urandom(16))' > instance/config.py`.  Defaults to the key `dev`, **should be explicitly set in a production environment**.
* `ALWAYS_REPORT_HTTPS` - If True, VExD will display addresses on the API page with the https scheme, regardless of what the web server thinks the correct scheme should be.  Helpful if there's a gateway in front of the web site.  Defaults to False
* `DOWNLOAD_DIRECTORY` - Where should VExD look for the statically created files on the Downloads page?  Defaults to `instance/downloads`
* `MONGODB_HOST`, `MONGODB_PORT` - The location of the mongo server, default to `localhost` and `27017`, respectively
* `MONGODB_USER`, `MONGODB_PASS` - The credentials for the mongo server.  If not provided, no credentials are sent to mongo

Any other [Flask configuration variable](https://flask.palletsprojects.com/en/2.0.x/config/?highlight=secret_key#builtin-configuration-values) can be set, though no others have been tested.

## Precomputing Results

Once VExD has been installed and configured (including the database) **and any time the database changes** be sure to (re-)generate the static database dumps available on the Downloads page, as well as precompute the background distribution used for statistics on the Enrichment page.  This can be done from the vexd source directory:

```
$ source env/bin/activate     # If necessary
(env) $ flask create-downloads
(env) $ flask save-distribution
```

## Running

### Development

In a development environment, it is easiest to [directly run a Flask server](https://flask.palletsprojects.com/en/2.0.x/server/) on http://localhost:5000.

```
$ source env/bin/activate                     # If necessary
(env) $ echo "FLASK_ENV = development" > .env # Only needs to be done once
(env) $ flask run
```

### Production

In a production environment, any WSGI container should work.  I have only used Apache and [mod_wsgi](https://modwsgi.readthedocs.io/en/master/).  Be sure to read the documentation for both of those programs thoroughly, as there are a number of options and security concerns!  However, the key steps are to:

1. Install VExD (following the above instructions) in some directory outside the web root.  I used `/srv/vexd`.
2. Create a `vexd.wsgi` file containing the following:
   ```
   from vexd import create_app

   application = create_app()
   ```
3. In your Apache config, include the following lines:
   ```
    WSGIDaemonProcess vexd python-home=/srv/vexd/env
    WSGIProcessGroup vexd
    WSGIApplicationGroup %{GLOBAL}
    WSGIScriptAlias / /var/www/wsgi/vexd.wsgi
   ```

The `WSGIScriptAlias` directive can also be used to serve VExD out of a subdirectory.  The `python-home` element of the `WSGIDaemonProcess` directive points to the virtual environment where VExD was installed.  With this configuration, you can update a running server by simply running `touch /var/www/wsgi/vexd.wsgi`, as opposed to reloading the entire Apache server.

These Apache directives are not sufficient by themselves (you still need to allow the server to access the `/var/www/wsgi` directory, for instance), but they are the core of the configuration.

## Javascript Libraries

VExD makes use of several open-source Javascript libraries, all of which are included in the `static` directory:

* [OLS autocomplete](https://www.npmjs.com/package/ols-autocomplete), from the [Ontology Lookup Search](https://www.ebi.ac.uk/ols/) at EBI.  Used for the BTO term lookup
* [OLS treeview](https://www.npmjs.com/package/ols-treeview), also from OLS.  Used to display the BTO hierarchy
* [Sortable](https://github.hubspot.com/sortable/).  Used to make the gene result table sortable

The two OLS libraries have quite a few dependencies ([JQuery](https://jquery.com), [Handlebars](https://handlebarsjs.com), the core-js fork of [typeahead.js](https://typeahead.js.org/), and [jsTree](https://www.jstree.com)).  These are all linked from [cdnjs](https://cdnjs.com) instead of being separately downloaded and included with VExD.