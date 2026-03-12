from django.conf import settings
from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ("core", "0002_sequencefeature"),
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    operations = [
        migrations.CreateModel(
            name="PCRProduct",
            fields=[
                ("id", models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name="ID")),
                ("name", models.CharField(max_length=255)),
                ("record_id", models.CharField(max_length=255)),
                ("forward_primer_label", models.CharField(blank=True, default="", max_length=255)),
                ("reverse_primer_label", models.CharField(blank=True, default="", max_length=255)),
                ("start", models.IntegerField()),
                ("end", models.IntegerField()),
                ("sequence", models.TextField()),
                ("length", models.IntegerField(default=0)),
                ("created_at", models.DateTimeField(auto_now_add=True)),
                ("creator", models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name="pcrproduct_created", to=settings.AUTH_USER_MODEL)),
                ("forward_feature", models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, related_name="pcr_products_as_forward_feature", to="core.sequencefeature")),
                ("forward_primer", models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, related_name="pcr_products_as_forward", to="core.primer")),
                ("reverse_feature", models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, related_name="pcr_products_as_reverse_feature", to="core.sequencefeature")),
                ("reverse_primer", models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, related_name="pcr_products_as_reverse", to="core.primer")),
                ("sequence_file", models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name="pcr_products", to="core.sequencefile")),
                ("users", models.ManyToManyField(blank=True, related_name="pcrproduct_access", to=settings.AUTH_USER_MODEL)),
            ],
            options={
                "ordering": ["-created_at", "name"],
            },
        ),
    ]
