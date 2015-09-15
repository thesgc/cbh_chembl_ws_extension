# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_core_model', '0017_auto_20150915_0022'),
    ]

    operations = [
        migrations.CreateModel(
            name='ChemregProject',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('cbh_core_model.project',),
        ),
    ]
