from django.conf import settings
from django.db import migrations, models
import django.db.models.deletion


def grant_sequencefile_access_to_uploaders(apps, schema_editor):
    SequenceFile = apps.get_model("core", "SequenceFile")
    for sequence_file in SequenceFile.objects.all():
        sequence_file.users.add(sequence_file.uploaded_by_id)


class Migration(migrations.Migration):
    dependencies = [
        ("core", "0004_project_pcr_products"),
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    operations = [
        migrations.AddField(
            model_name="sequencefile",
            name="users",
            field=models.ManyToManyField(blank=True, related_name="sequencefile_access", to=settings.AUTH_USER_MODEL),
        ),
        migrations.CreateModel(
            name="AnalysisJob",
            fields=[
                ("id", models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name="ID")),
                (
                    "job_type",
                    models.CharField(
                        choices=[
                            ("primer_binding", "Primer binding"),
                            ("pcr_product_discovery", "PCR product discovery"),
                        ],
                        max_length=50,
                    ),
                ),
                ("celery_task_id", models.CharField(blank=True, default="", max_length=255)),
                (
                    "status",
                    models.CharField(
                        choices=[
                            ("pending", "Pending"),
                            ("running", "Running"),
                            ("success", "Success"),
                            ("failure", "Failure"),
                        ],
                        default="pending",
                        max_length=20,
                    ),
                ),
                ("error_message", models.TextField(blank=True, default="")),
                ("result_payload", models.JSONField(blank=True, null=True)),
                ("target_object_id", models.PositiveBigIntegerField(blank=True, null=True)),
                ("created_at", models.DateTimeField(auto_now_add=True)),
                ("updated_at", models.DateTimeField(auto_now=True)),
                (
                    "owner",
                    models.ForeignKey(
                        on_delete=django.db.models.deletion.CASCADE,
                        related_name="analysis_jobs",
                        to=settings.AUTH_USER_MODEL,
                    ),
                ),
                (
                    "primer",
                    models.ForeignKey(
                        blank=True,
                        null=True,
                        on_delete=django.db.models.deletion.CASCADE,
                        related_name="analysis_jobs",
                        to="core.primer",
                    ),
                ),
                (
                    "primer_pair",
                    models.ForeignKey(
                        blank=True,
                        null=True,
                        on_delete=django.db.models.deletion.CASCADE,
                        related_name="analysis_jobs",
                        to="core.primerpair",
                    ),
                ),
                (
                    "sequence_file",
                    models.ForeignKey(
                        blank=True,
                        null=True,
                        on_delete=django.db.models.deletion.CASCADE,
                        related_name="analysis_jobs",
                        to="core.sequencefile",
                    ),
                ),
            ],
            options={"ordering": ["-created_at"]},
        ),
        migrations.RunPython(
            grant_sequencefile_access_to_uploaders,
            migrations.RunPython.noop,
        ),
    ]
