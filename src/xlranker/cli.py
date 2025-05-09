import cyclopts
import xlranker.ml.data as xlr_data

app = cyclopts.App()


@app.command()
def hello(name: str):
    """greets the user!"""
    print(f"Hello {name}!")


@app.command()
def test_loading():
    df = xlr_data.load_default_ppi()
    print(df.head())


@app.command()
def start(input: str):
    """Run the full prioritization pipeline

    Requires input file to be in the format specified in the project documentation.

    Args:
        input (str): input file in TSV format
    """
    pass


def cli():
    app()
