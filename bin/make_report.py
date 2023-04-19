import json
import csv
import collections
from mako.lookup import TemplateLookup
import datetime as dt
from datetime import date

from mako.template import Template
from mako.runtime import Context
from mako.exceptions import RichTraceback
from io import StringIO
import os
from Bio import SeqIO


def make_output_report(report_to_generate):

    # collate data for tables in the report
    data_for_report = {"KEY_SUMMARY_TABLE":[],"KEY_COMPOSITION_TABLE":[]}
    print("have data")

    template_dir = os.path.abspath(os.path.dirname('bin'))

    mylookup = TemplateLookup(directories=["bin"]) #absolute or relative works
    print(mylookup)
    print("found template")
    mytemplate = mylookup.get_template("pantheon_report.mako")
    print(mytemplate)
    buf = StringIO()

    ctx = Context(buf,
                    date = date.today(),
                    run_name="test_runname",
                    version = "__version__",
                    sample="sample",
                    data_for_report = data_for_report)
    print("context")

    try:
        mytemplate.render_context(ctx)
        print("render")

    except:
        print("failed")
        traceback = RichTraceback()
        for (filename, lineno, function, line) in traceback.traceback:
            print("File %s, line %s, in %s" % (filename, lineno, function))
            print(line, "\n")
        print("%s: %s" % (str(traceback.error.__class__.__name__), traceback.error))

    with open(report_to_generate, 'w') as fw:
        print("Generating: " + f"{report_to_generate}")
        fw.write(buf.getvalue())


def main():
    print("starting")
    make_output_report("out.html")


if __name__ == "__main__":
    main()