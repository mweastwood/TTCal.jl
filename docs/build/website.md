
<a id='Building-This-Website-1'></a>

# Building This Website


This website is built using [MkDocs](http://www.mkdocs.org/) and hosted on [Github Pages](https://pages.github.com/). Building and deploying the documentation will require a little bit of setup.


1. Run `pip install --user mkdocs` to obtain MkDocs.
2. Run `pip install --user python-markdown-math` to obtain the markdown    extension that allows LaTeX-style math.
3. Edit the documentation in the `docs/` directory.
4. Run `mkdocs serve` and navigate your web browser to    [http://127.0.0.1:8000](http://127.0.0.1:8000)    to check that everything looks correct.
5. Commit all your changes.
6. Run `mkdocs gh-deploy` to deploy the new documentation to Github Pages.

