# Publishing a PyMIRAGER 0.x release

The source bundle is prepared offline; it is **not itself a commit, tag or GitHub release**.

1. Inspect your repository (`git status`, `git branch --show-current`, `git remote -v`) and review the update bundle before copying it. If the remote repository has newer changes, merge them deliberately rather than replacing them blindly.
2. Copy the contents of the v0.4.7 update archive into the **root of your GitHub checkout**, preserving paths.
3. Run the local tests and verify metadata:

   ```bash
   python3 -m pip install -e .
   python3 -m unittest discover -s tests -v
   python3 -m jupiter_rt --version
   git diff --check
   git diff
   ```

4. Check the Jupiter input case and external opacities; rerun the grey-aerosol case using `--preset aerosol_cutoff`, which now allows **85–119**.
5. Commit, tag and push from your own authenticated GitHub checkout:

   ```bash
   git add README.md CHANGELOG.md pyproject.toml jupiter_rt/cli.py jupiter_rt/__init__.py docs/ tests/
   git commit -m "Release PyMIRAGER v0.4.7"
   git tag -a v0.4.7 -m "PyMIRAGER v0.4.7"
   git push origin HEAD
   git push origin v0.4.7
   ```

6. Optionally create a GitHub Release from tag `v0.4.7`, using `docs/RELEASE_NOTES_0.4.7.md` as its description.

**Do not tag 0.4.7 before checking that the remote matches the intended source, because the release archive was generated from the v0.4.6 development snapshot rather than by authenticated access to the current GitHub HEAD.**
