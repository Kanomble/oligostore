def assign_creator(obj, user):
    obj.creator = user
    obj.save()
    obj.users.add(user)
    return obj
