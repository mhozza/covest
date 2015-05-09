#! /usr/bin/env python
import pystache


def lines_to_list(header, lines):
    res = list()
    for l in lines:
        res.append(
            [l[k] for k in header]
        )
    return res


def format_table(header, lines, template_file, escape=None, round_floats=4):
    def format(val):
        if round_floats and type(val) is float:
            val = round(val, round_floats)
        if escape:
            try:
                return escape(val)
            except:
                pass
        return val
    lines = lines_to_list(header, lines)
    data = {
        'header': [
            {'value': format(h), 'first': i == 0, 'last': i == len(header) - 1}
            for i, h in enumerate(header)
        ],
        'body': [
            {'line': [
                {'value': format(v), 'first': i == 0, 'last': i == len(l) - 1}
                for i, v in enumerate(l)
            ]} for l in lines
        ],
    }

    with open(template_file) as f:
        template = f.read()
        return pystache.render(template, data)
