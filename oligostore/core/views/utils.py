from django.core.paginator import Paginator


def paginate_queryset(request, qs, per_page=10):
    paginator = Paginator(qs, per_page)
    page_number = request.GET.get("page")
    page_obj = paginator.get_page(page_number)
    query_params = request.GET.copy()
    query_params.pop("page", None)
    return page_obj, query_params.urlencode()