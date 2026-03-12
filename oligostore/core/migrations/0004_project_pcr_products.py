from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ("core", "0003_pcrproduct"),
    ]

    operations = [
        migrations.AddField(
            model_name="project",
            name="pcr_products",
            field=models.ManyToManyField(blank=True, related_name="projects", to="core.pcrproduct"),
        ),
    ]
