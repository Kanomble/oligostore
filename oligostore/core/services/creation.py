from django.db import transaction

from ..models import PCRProduct, Primer, PrimerPair, SequenceFeature, SequenceFile
from .ownership import assign_creator, grant_user_access


@transaction.atomic
def create_owned_project(*, form, user):
    project = form.save(commit=False)
    assign_creator(project, user)
    project.save()
    grant_user_access(project, user)
    return project


@transaction.atomic
def create_owned_primer_pair(*, name, forward_primer, reverse_primer, user):
    pair = PrimerPair(
        name=name,
        forward_primer=forward_primer,
        reverse_primer=reverse_primer,
    )
    assign_creator(pair, user)
    pair.save()
    grant_user_access(pair, user)
    return pair


@transaction.atomic
def create_primer_pair_with_new_primers(
    *,
    pair_name,
    forward_name,
    forward_sequence,
    reverse_name,
    reverse_sequence,
    user,
    forward_overhang="",
    reverse_overhang="",
):
    forward = Primer.create_with_analysis(
        primer_name=forward_name,
        sequence=forward_sequence,
        overhang_sequence=forward_overhang,
        user=user,
    )
    reverse = Primer.create_with_analysis(
        primer_name=reverse_name,
        sequence=reverse_sequence,
        overhang_sequence=reverse_overhang,
        user=user,
    )
    pair = create_owned_primer_pair(
        name=pair_name,
        forward_primer=forward,
        reverse_primer=reverse,
        user=user,
    )
    return forward, reverse, pair


@transaction.atomic
def create_primer_and_optional_feature(
    *,
    user,
    sequence_file: SequenceFile,
    primer_name: str,
    sequence: str,
    overhang_sequence: str,
    save_to_primers: bool,
    feature_attachment=None,
):
    primer = None
    if save_to_primers:
        primer = Primer.create_with_analysis(
            primer_name=primer_name,
            sequence=sequence,
            overhang_sequence=overhang_sequence,
            user=user,
        )

    attached_feature = None
    if feature_attachment:
        attached_feature = SequenceFeature.objects.create(
            sequence_file=sequence_file,
            primer=primer,
            record_id=feature_attachment["record_id"],
            start=feature_attachment["start"],
            end=feature_attachment["end"],
            strand=feature_attachment["strand"],
            feature_type=SequenceFeature.TYPE_PRIMER_BIND,
            label=primer.primer_name if primer else primer_name,
            created_by=user,
        )

    return primer, attached_feature


@transaction.atomic
def create_pcr_product(
    *,
    user,
    sequence_file,
    name,
    record_id,
    start,
    end,
    sequence,
    forward_primer=None,
    reverse_primer=None,
    forward_feature=None,
    reverse_feature=None,
    forward_primer_label="",
    reverse_primer_label="",
):
    product = PCRProduct(
        name=name,
        sequence_file=sequence_file,
        record_id=record_id,
        forward_primer=forward_primer,
        reverse_primer=reverse_primer,
        forward_feature=forward_feature,
        reverse_feature=reverse_feature,
        forward_primer_label=forward_primer_label or getattr(forward_primer, "primer_name", ""),
        reverse_primer_label=reverse_primer_label or getattr(reverse_primer, "primer_name", ""),
        start=start,
        end=end,
        sequence=sequence,
    )
    assign_creator(product, user)
    product.save()
    grant_user_access(product, user)
    return product
