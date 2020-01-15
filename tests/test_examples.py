def test_examples():
    from os import listdir
    from os.path import isfile, join
    import subprocess

    path = "examples/"

    example_files = [f for f in listdir(path) if isfile(join(path, f))]
    example_files = [f for f in example_files if f[-3:] == ".py"]

    for example in example_files:
        return_code = subprocess.call(['python', path+example])
        assert not return_code


if __name__ == "__main__":
    test_examples()
