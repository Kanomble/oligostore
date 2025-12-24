from django import template

register = template.Library()

@register.filter(name='add_class')
def add_class(field, css_class):
    """Adds a CSS class to a Django form field widget."""
    return field.as_widget(attrs={
        **field.field.widget.attrs,
        "class": css_class
    })

@register.filter(name='get_item')
def get_item(obj, key):
    """Allows dynamic lookup: form|get_item:'FIELD_NAME'."""
    return obj[key]