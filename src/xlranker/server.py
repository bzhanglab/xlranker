import importlib.resources as pkg_resources
import mimetypes
import threading
import time
import uuid
import webbrowser
from importlib.resources import files

from flask import Flask, Response, jsonify, request

import xlranker.static
from xlranker.util import map_col_sep
from xlranker.util.readers import read_network

app = Flask(__name__, static_folder=str(files("xlranker.static")))


# Store tasks in memory (for demo; in production use Redis/DB)
tasks = {}


def long_task(task_id, name):
    # simulate long 5-min task
    time.sleep(10)
    tasks[task_id]["status"] = "done"
    tasks[task_id]["result"] = f"Task for {name} finished!"


@app.route("/")
def index():
    data = pkg_resources.files(xlranker.static).joinpath("index.html").read_bytes()
    return Response(data, mimetype="text/html")


@app.route("/<path:filename>")
def serve_static(filename):
    path = pkg_resources.files(xlranker.static).joinpath(filename)
    if not path.is_file():
        return "Not Found", 404
    data = path.read_bytes()
    mimetype, _ = mimetypes.guess_type(filename)
    return Response(data, mimetype=mimetype or "application/octet-stream")


@app.route("/api/upload_network", methods=["POST"])
def upload_network():
    try:
        uploaded_file = request.files.get("peptide_network")
        req = request.form.to_dict()
        sep = map_col_sep(req["col_sep"])
        if not uploaded_file:
            return jsonify({"status": "err", "message": "No file uploaded"}), 400
        res = read_network(uploaded_file.stream, col_sep=sep)
        if res.duplicate_rows > 0:  # Send warning that duplicate edges were removed.
            return jsonify(
                {
                    "status": "warn",
                    "message": f"Found and removed {res.duplicate_rows} duplicated edge(s) in network.",
                    "filename": uploaded_file.filename,
                    "sep": sep,
                    "duplicates": res.duplicate_rows,
                }
            )
        return jsonify(
            {
                "status": "ok",
                "message": f"File processed successfully with {len(res.network)} edges.",
            }
        )
    except Exception as e:
        return jsonify({"status": "err", "message": str(e)}), 500


@app.route("/start_task", methods=["POST"])
def start_task():
    data = request.get_json()
    name = data.get("name", "guest")

    task_id = str(uuid.uuid4())
    tasks[task_id] = {"status": "running", "result": None}

    # Run in background thread
    thread = threading.Thread(target=long_task, args=(task_id, name))
    thread.start()

    return jsonify({"task_id": task_id})


@app.route("/task/<task_id>", methods=["GET"])
def get_task(task_id):
    task = tasks.get(task_id)
    if not task:
        return jsonify({"error": "Task not found"}), 404
    return jsonify(task)


def start_server(host: str = "127.0.0.1", port: str = "5000", debug: bool = False):
    url = f"http://{host}:{port}/"
    # open browser in another thread so it doesn't block Flask
    threading.Timer(1.0, lambda: webbrowser.open(url)).start()
    app.run(debug=debug)
