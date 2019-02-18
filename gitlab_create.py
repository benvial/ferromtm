#!/usr/bin/env python

import gitlab
gl = gitlab.Gitlab.from_config('gitlab')
group_id = 4428745
project = gl.projects.create({'name': 'ferromtm', 'namespace_id': group_id})
