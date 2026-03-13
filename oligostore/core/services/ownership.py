from django.db import models


def assign_creator(obj: models.Model, user):
    obj.creator = user
    return obj


def grant_user_access(obj: models.Model, user):
    obj.users.add(user)
    return obj
