from behave import given, when, then
import json
from django.contrib.auth.models import User, Group

@given('I have a valid molfile')
def step(context):
    context.post_data["ctab"] = """


  8  8  0  0  0  0            999 V2000
    0.0000    1.0000    0.0000 C   0  0  0  0  0  0
    0.8660    0.5000    0.0000 C   0  0  0  0  0  0
    0.8660   -0.5000    0.0000 C   0  0  0  0  0  0
    0.0000   -1.0000    0.0000 C   0  0  0  0  0  0
   -0.8660   -0.5000    0.0000 C   0  0  0  0  0  0
   -0.8660    0.5000    0.0000 C   0  0  0  0  0  0
    0.0000    2.0000    0.0000 O   0  0  0  0  0  0
    0.0000   -2.0000    0.0000 O   0  0  0  0  0  0
  1  21  0     0  0
  2  32  0     0  0
  3  41  0     0  0
  4  51  0     0  0
  5  62  0     0  0
  6  11  0     0  0
  1  72  0     0  0
  4  82  0     0  0
M  END"""

@given("my user is member of a group")
def step(context):
    context.g = Group.objects.create(name="g")
    context.g.user_set.add(context.u)

@given("a valid project exists {projkey}")
def step(context, projkey):
    from cbh_chembl_model_extension.models import Project
    Project.objects.create(name=projkey, project_key=projkey, created_by=context.u)

@given("I automatically have editor permissions as creator")
def step(context):
    from cbh_chembl_model_extension.models import Project
    for p in Project.objects.all():
        assert p.created_by.has_perm("%d.%s" % ( p.id, "editor"))

@given("I remove my permissions")
def step(context):
    for perm in context.u.user_permissions.all():
        context.u.user_permissions.remove(perm)
        context.u = User.objects.get(pk=context.u.id)
    from cbh_chembl_model_extension.models import Project
    for p in Project.objects.all():
        assert p.created_by.has_perm("%d.%s" % ( p.id, "editor")) == False


@given("I have given {me_or_group} editor privileges for {projkey}")
def step(context,me_or_group, projkey):
    from cbh_chembl_model_extension.models import Project
    p = Project.objects.get(project_key=projkey)
    if me_or_group == "myself":
        p.make_editor(context.u)
    elif me_or_group == "mygroup":
        p.make_editor(context.g)
    context.u = User.objects.get(pk=context.u.id)
    assert context.u.has_perm("%d.%s" % ( p.id, "editor"))


@given("I have given {me_or_group} viewer privileges for {projkey}")
def step(context,me_or_group, projkey):
    from cbh_chembl_model_extension.models import Project
    p = Project.objects.get(project_key=projkey)
    if me_or_group == "myself":
        p.make_viewer(context.u)
    elif me_or_group == "mygroup":
        p.make_viewer(context.g)
    context.u = User.objects.get(pk=context.u.id)
    assert context.u.has_perm("%d.%s" % ( p.id, "viewer"))


@given("I havent given {grouporuser} {editororviewer} privileges for {projkey}")
def step(context,grouporuser, editororviewer, projkey):
    pass


@given("I have a valid list of SMILES")
def step(context):
    pass


@given("I have an invalid Excel File")
def step(context):
    pass


@given("a single batch exists in {projkey}")
def step(context, projkey):
    from cbh_chembl_model_extension.models import Project, CBHCompoundBatch
    p = Project.objects.get(project_key=projkey)
    context.batch = CBHCompoundBatch.objects.create(ctab=context.post_data["ctab"], project=p)

@then("I {action} my cbh_compound_batch to {projkey} and the response code will be {responsecode}")
def step(context, action=None, projkey=None, responsecode=None):
    from cbh_chembl_model_extension.models import Project    
    if action in ["validate","create"]:
        if action =="validate":
            path = "/chemblws/cbh_compound_batches/validate/"
            func = context.api_client.post
        elif action == "create":
            path = "/chemblws/cbh_compound_batches/"
            func = context.api_client.post
        context.post_data["projectKey"] = projkey
        resp = func(
            path,
            format='json' ,
            data=context.post_data,
            )
        assert resp.status_code == int(responsecode)
    else:
        from cbh_chembl_model_extension.models import Project
        
        path = "/chemblws/cbh_compound_batches/"
        if action == "get":
            path = path + str(context.batch.id)
            print path
        resp = context.api_client.get(
            path,
            )
        if action == "get":
            assert resp.status_code == int(responsecode)
        else:
            project_id = Project.objects.get(project_key=projkey).id
           
            assert len(context.ser.deserialize(resp.content)["objects"]) == int(responsecode.split("_")[0])





