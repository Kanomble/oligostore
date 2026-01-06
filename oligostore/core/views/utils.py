from django.core.paginator import Paginator


def paginate_queryset(request, qs, per_page=10):
    allowed_per_page = {10, 20, 50, 100}
    requested_per_page = request.GET.get("per_page")
    if requested_per_page:
        try:
            requested_per_page = int(requested_per_page)
        except (TypeError, ValueError):
            requested_per_page = None
    if requested_per_page in allowed_per_page:
        per_page = requested_per_page

    paginator = Paginator(qs, per_page)
    page_number = request.GET.get("page")
    page_obj = paginator.get_page(page_number)
    query_params = request.GET.copy()
    query_params.pop("page", None)
    return page_obj, query_params.urlencode()