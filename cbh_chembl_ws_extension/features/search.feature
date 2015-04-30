
Feature: Search results are consistent with the search parameters
'''
    I wish to search for a compound or set of compounds constrained by various parameters

    Scenario: User submits a search query with no constraints.
        Given a User
        and my user is a member of a group
        and I have viewer permissions
        When I log in
        #search parameter building statements would go at this point
        then I submit my search
        and I can see search results
        and all the results are viewable by user

    Scenario: User submits a search query with date constraints.
        Given a User
        and my user is a member of a group
        and I have given myself viewer privileges for proja
        When I log in
        and I specify a start date start
        and I specify an end date end
        then I submit my search
        and I can see search results
        and all the results are viewable by user
        and all the results were submitted after start
        and all the results were submitted before end

    Scenario: User submits a search query with a single project constraint.
        Given a User
        and my user is a member of a group
        and I have given myself viewer privileges for proja
        and I have viewer permissions
        When I log in
        and I specify a project proja
        then I submit my search
        and I can see search results
        and all the results are viewable by user
        and all the results belong to proja

    Scenario: User submits a search query with a custom field constraint.
        Given a User
        and my user is a member of a group
        and I have given myself viewer privileges for proja
        and I have viewer permissions
        When I log in
        and I specify a project proja
        and I specify a custom field tag taga
        then I submit my search
        and I can see search results
        and all the results are viewable by user
        and all the results belong to proja
        and all the results contain taga

    Scenario: User submits a search query with a structure to find a structural exact match.
        Given a User
        and my user is a member of a group
        and I have given myself viewer privileges for proja
        and I have viewer permissions
        When I log in
        and I specify a project proja
        and I specify a structure struca
        and I specify structural search type exact
        then I submit my search
        and I can see search results
        and all the results are viewable by user
        and all the results belong to proja
        and all the results exact match struca

    Scenario: User submits a search query with a structure to find a structural substructure match.
        Given a User
        and my user is a member of a group
        and I have given myself viewer privileges for proja
        and I have viewer permissions
        When I log in
        and I specify a project proja
        and I specify a structure struca
        and I specify structural search type substructure
        then I submit my search
        and I can see search results
        and all the results are viewable by user
        and all the results belong to proja
        and all the results substructure match struca
