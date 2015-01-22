Feature: CBH Compound Batch
    In order to keep track of compounds I would like to create 
    a batch of a new compound in the project I am an admin of.
    As a project stakeholder I do not want batches to be logged in 
    projects that the user does not have editor access to


    Scenario Outline: Batches project privileges
        Given a User
        When <dologin> log in
        and I have a valid molfile
        and I have a valid project key proja
        and I have <editor> editor privileges for proja
        and I have <viewer> viewer privileges for proja
        when I save my cbh_compound_batch to proja
        Then the response code will be <responsecode>

        Examples: Tests
        |   dologin    |   editor  |   viewer  |   responsecode   | 
        |      I         |           |           |   201             |
        |       I        |           |   not     |   201             |
        |        I       |   not     |           |   401             |
        |         I      |   not     |   not     |   401             |
        |   do not      |           |           |   401             |
 

    """
    Created using
    arrays = [("","do not"), ("","not"), ("","not"), ("401",)]
    rows = list(itertools.product(*arrays))
    rows = ["|\t" +  "\t\t|\t".join(row) + "\t\t|" for row in   list(itertools.product(*arrays))]
    print "\n".join(rows)
    """





    #Scenario Outline: Other project privileges, ivladifd project key

