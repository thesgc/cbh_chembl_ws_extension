# -*- coding: utf-8 -*-
"""steps/browser_steps.py -- step implementation for our browser feature demonstration.
"""
from behave import given, when, then
import json

 
@given('a user')
def step(context):
    from django.contrib.auth.models import User
    u = User(username='foo', email='foo@example.com')
    u.set_password('bar')
    u.save()
    context.u = u
 


@when('I do not log in')
def step(context):
    pass
 
 

 
@when('I log in')
def step(context):
    context.api_client.client.login(username="foo", password="bar")
 
 
@then('I see a 401 error')
def step(context):

    resp = context.api_client.get("/dev/users/",  format='json')
    
    assert resp.status_code == 401
