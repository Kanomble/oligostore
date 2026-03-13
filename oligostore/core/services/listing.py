from django.db.models import Q


def apply_search(queryset, query, fields):
    if not query:
        return queryset

    clause = Q()
    for field in fields:
        clause |= Q(**{f"{field}__icontains": query})
    return queryset.filter(clause)


def apply_ordering(queryset, requested_order, allowed_orders, default_order):
    return queryset.order_by(allowed_orders.get(requested_order, default_order))
