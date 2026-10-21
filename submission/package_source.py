"""Build the final PDF and source ZIP without changing manuscript or figures."""

from collections import Counter
from pathlib import Path
import re
import shutil
import subprocess
from tempfile import TemporaryDirectory
from zipfile import ZipFile, ZIP_DEFLATED


ROOT = Path(__file__).resolve().parents[1]
MAIN = ROOT / "doc" / "main"
ARCHIVE = ROOT / "submission" / "jgaa-source.zip"
MANUSCRIPT = MAIN / "jgaa_without_red.tex"


def source_files() -> set[Path]:
    files = {MANUSCRIPT}
    files.update(MAIN / name for name in (
        "jgaa-art.cls", "abbrvurl.bst", "by.eps", "by.png",
        "orcid.eps", "orcid.png", "cube-spread.pdf", "cube-spread.ps",
        "FruchtermanReingoldByRandomSubspace.bib",
        "individual/individual.ipynb", "individual/legend.py",
    ))
    files.update((MAIN / "individual" / "hist").glob("*.txt"))
    pending = [MANUSCRIPT]
    visited = set()
    while pending:
        tex = pending.pop()
        if tex in visited:
            continue
        visited.add(tex)
        text = re.sub(r"(?m)(?<!\\)%.*$", "", tex.read_text(encoding="utf-8"))
        for name in re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", text):
            graphic = tex.parent / name
            if not graphic.is_file():
                graphic = MAIN / name
            if not graphic.is_file():
                raise FileNotFoundError(graphic)
            graphic = graphic.resolve()
            files.add(graphic)
            for extension in (".tex", ".py", ".ipynb"):
                sidecar = graphic.with_suffix(extension)
                if sidecar.is_file():
                    files.add(sidecar)
                    if extension == ".tex":
                        pending.append(sidecar)
    for path in files:
        path.resolve().relative_to(MAIN.resolve())
        if not path.is_file():
            raise FileNotFoundError(path)
    return files


def main() -> None:
    bib = (MAIN / "FruchtermanReingoldByRandomSubspace.bib").read_text(encoding="utf-8")
    keys = re.findall(r"(?m)^@\w+\s*\{\s*([^,\s]+)\s*,", bib)
    duplicates = [key for key, count in Counter(keys).items() if count > 1]
    if duplicates:
        raise ValueError(f"Duplicate bibliography keys: {', '.join(duplicates)}")

    files = source_files()
    # JGAA requests complete LaTeX sources: https://jgaa.info/index.php/jgaa/about/submissions
    with TemporaryDirectory(prefix="jgaa-source-") as directory:
        temporary = Path(directory)
        build = temporary / "jgaa-source"
        for source in sorted(files):
            name = Path("jgaa.tex") if source == MANUSCRIPT else source.relative_to(MAIN)
            target = build / name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, target)
        result = subprocess.run(
            ["latexmk", "-norc", "-pdf", "-interaction=nonstopmode",
             "-halt-on-error", "-file-line-error", "jgaa.tex"],
            cwd=build, capture_output=True, text=True, errors="replace",
        )
        if result.returncode:
            print(result.stdout[-12000:])
            print(result.stderr)
            result.check_returncode()
        log = (build / "jgaa.log").read_text(encoding="utf-8", errors="replace")
        for warning in ("Overfull", "Float too large", "undefined references", "undefined citations"):
            if warning in log:
                raise RuntimeError(f"Check compilation warning: {warning}")
        archive = temporary / "jgaa-source.zip"
        with ZipFile(archive, "w", compression=ZIP_DEFLATED) as output:
            for source in sorted(files):
                name = Path("jgaa.tex") if source == MANUSCRIPT else source.relative_to(MAIN)
                output.write(build / name, "jgaa-source/" + name.as_posix())
        with ZipFile(archive) as output:
            if output.testzip() is not None:
                raise RuntimeError("ZIP integrity check failed")
            for source in files:
                name = Path("jgaa.tex") if source == MANUSCRIPT else source.relative_to(MAIN)
                if output.read("jgaa-source/" + name.as_posix()) != source.read_bytes():
                    raise RuntimeError(f"ZIP content differs: {source}")
        shutil.copyfile(build / "jgaa.pdf", MANUSCRIPT.with_suffix(".pdf"))
        shutil.copyfile(archive, ARCHIVE)
    print(f"Built PDF and {ARCHIVE} ({len(files)} source files).")


if __name__ == "__main__":
    main()
