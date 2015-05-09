#! /usr/bin/env python
from collections import defaultdict


def combine_tables(lines1, lines2, template='{} \u00B1 {}', round_floats=4):
    lines = list()
    for l1, l2 in zip(lines1, lines2):
        l = defaultdict(str)
        for key, val1 in l1.items():
            val2 = l2[key]
            if round_floats:
                if type(val1) is float:
                    val1 = round(val1, round_floats)
                if type(val2) is float:
                    val2 = round(val2, round_floats)
            l[key] = template.format(val1, val2)
        lines.append(l)
    return lines


def format_table(header, lines, template, line_template, header_cell_template,
                 cell_template, round_floats=4):
    lbody = ''
    for c in header:
        lbody += header_cell_template.format(cbody=c)
    thead = line_template.format(lbody=lbody)
    tbody = ''
    for l in lines:
        lbody = ''
        for c in header:
            val = l.get(c, None)
            if round_floats and type(val) is float:
                val = round(val, round_floats)
            lbody += cell_template.format(cbody=val)
        tbody += line_template.format(lbody=lbody)
    return template.format(thead=thead, tbody=tbody)


def format_table_html(header, lines):
    template = '''<!DOCTYPE html>
<html>
<head>
<!-- Latest compiled and minified CSS -->
<link rel="stylesheet"
  href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.4/css/bootstrap.min.css">
</head>
<table class="table table-condensed table-hover">
  {thead}
  {tbody}
</table>
</html>
    '''
    line_template = '<tr>\n{lbody}  </tr>'
    cell_template = '    <td>{cbody}</td>\n'
    header_cell_template = '    <th>{cbody}</th>\n'
    return format_table(header, lines, template, line_template,
                        header_cell_template, cell_template)


def format_table_csv(header, lines):
    template = '{thead}\n{tbody}'
    line_template = '{lbody}\n'
    cell_template = '{cbody}, '
    header_cell_template = cell_template
    return format_table(header, lines, template, line_template,
                        header_cell_template, cell_template)
