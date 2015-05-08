@then("I get projects and the response code will be {responsecode} with {project} in it")
def step(context, responsecode, project):

    resp = context.api_client.get("/dev/cbh_projects/",  format='json')
    print(resp.status_code)
    assert int(resp.status_code) == int(responsecode)
    assert project in resp.content


@then("I get projects and the response code will be {responsecode} without {project} in it")
def step(context, responsecode, project):

    resp = context.api_client.get("/dev/cbh_projects/",  format='json')
    print(resp.status_code)
    assert int(resp.status_code) == int(responsecode)
    assert project not in resp.content


@then("I get projects and the response code will be {responsecode}")
def step(context, responsecode):
    resp = context.api_client.get("/dev/cbh_projects/",  format='json')
    assert int(resp.status_code) == int(responsecode)
