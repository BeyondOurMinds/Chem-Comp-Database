def test_cli_imports():
    # Basic smoke test to ensure the CLI entry point imports cleanly.
    from Chem_Database.gui import app

    assert app is not None
