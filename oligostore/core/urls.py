from django.urls import path
from .views import primer_list, primer_create,\
    primer_delete, home, register, primerpair_list, \
    primerpair_create, project_list,\
    primerpair_combined_create, primerpair_delete, \
    project_create, project_dashboard,\
    project_add_primerpair, project_remove_primerpair,\
    analyze_primer_view, analyze_primerpair_view,\
    analyze_sequence_view, save_generated_primerpair, \
    download_product_sequence, primer_binding_analysis, \
    sequencefile_list, sequencefile_upload, \
    project_primer_list, project_download_sequence_files, \
    primer_binding_analysis_async, primer_binding_status, \
    project_add_sequencefile, project_remove_sequencefile, \
    download_selected_primers, download_selected_primerpairs, \
    primer_import_excel

urlpatterns = [
    # Projects
    path("projects/", project_list, name="project_list"),
    path("projects/create/", project_create, name="project_create"),
    path("projects/<int:project_id>/", project_dashboard, name="project_dashboard"),
    path(
        "projects/<int:project_id>/primers/",
        project_primer_list,
        name="project_primer_list",
    ),
    path(
        "projects/<int:project_id>/sequence-files/download/",
        project_download_sequence_files,
        name="project_download_sequence_files",
    ),

    # Connecting primerpairs
    path("projects/<int:project_id>/add-pair/<int:pair_id>/", project_add_primerpair, name="project_add_primerpair"),
    path("projects/<int:project_id>/remove-pair/<int:pair_id>/", project_remove_primerpair, name="project_remove_primerpair"),
        path(
            "projects/<int:project_id>/add-sequence-file/<int:sequencefile_id>/",
            project_add_sequencefile,
            name="project_add_sequencefile",
        ),
        path(
            "projects/<int:project_id>/remove-sequence-file/<int:sequencefile_id>/",
            project_remove_sequencefile,
            name="project_remove_sequencefile",
        ),
    # Basic project paths
    path("", home, name="home"),
    path("register/", register, name="register"),
    path("primerpairs/delete/<int:primerpair_id>/", primerpair_delete, name="primerpair_delete"),
    path("primerpairs/create/combined/", primerpair_combined_create, name="primerpair_combined_create"),
    path("primerpair_create/", primerpair_create, name="primerpair_create"),
    path("primerpair_list",primerpair_list,name="primerpair_list"),
    path("project_list",project_list,name="project_list"),
    path("primer_list", primer_list, name="primer_list"),
    path("primer_list/import/", primer_import_excel, name="primer_import_excel"),
    path(
        "primer_list/download/",
        download_selected_primers,
        name="download_selected_primers",
    ),
    path("create/", primer_create, name="primer_create"),
    path(
        "primerpair_list/download/",
        download_selected_primerpairs,
        name="download_selected_primerpairs",
    ),
    path("delete/<int:primer_id>/", primer_delete, name="primer_delete"),

    # Sequence analysis
    path("analyze_sequence/", analyze_sequence_view, name="analyze_sequence"),
    path("save_primerpair/", save_generated_primerpair, name="save_generated_primerpair"),
    path(
        "primer-product/download/",
        download_product_sequence,
        name="download_product_sequence",
    ),
    path(
        "sequence-files/",
        sequencefile_list,
        name="sequencefile_list",
    ),
    path(
        "sequence-files/upload/",
        sequencefile_upload,
        name="sequencefile_upload",
    ),
    # Ajax primer analysis
    path("analyze-primer/", analyze_primer_view, name="analyze_primer"),
    path("analyze-primerpair/", analyze_primerpair_view, name="analyze_primerpair"),

    # Primer binding analysis
    path("primer-binding/", primer_binding_analysis, name="primer_binding"),
    path(
        "primer-binding/async/",
        primer_binding_analysis_async,
        name="primer_binding_async",
    ),
    path(
        "primer-binding/status/<str:task_id>/",
        primer_binding_status,
        name="primer_binding_status",
    ),
]
