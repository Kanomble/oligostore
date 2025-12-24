def assign_creator(obj, user):
    obj.creator = user  # 2
    obj.save()  # 3
    return obj
