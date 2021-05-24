"""

This is a performance test suitable for running with Scalene https://github.com/emeryberger/scalene

It is *not* run as a unit test

"""
import os

import puretabix

if __name__ == "__main__":
    targets = (
        ("1", 3538692),
        ("1", 6447088),
        ("1", 20916748),
        ("1", 25445572),
        ("1", 89619960),
        ("1", 109972419),
        ("1", 112960987),
        ("1", 113054364),
        ("1", 114807609),
        ("1", 148737614),
        ("1", 156527639),
        ("1", 202667513),
        ("1", 202691651),
        ("1", 208030764),
        ("1", 224113597),
        ("1", 233411738),
        ("2", 27754171),
        ("2", 154960831),
        ("2", 191960323),
        ("2", 219310743),
        ("2", 220087354),
        ("2", 230785920),
        ("2", 237952052),
        ("10", 50494625),
        ("10", 50497918),
        ("10", 50526658),
        ("10", 97164259),
        ("10", 114195064),
        ("10", 134015197),
        ("11", 5522482),
        ("11", 55814773),
        ("11", 87548096),
        ("11", 95217918),
        ("11", 107886290),
        ("11", 107888003),
        ("11", 111113468),
        ("12", 11066354),
        ("12", 25133862),
        ("12", 47455065),
        ("12", 93166154),
        ("12", 93183051),
        ("12", 107448080),
        ("12", 109821275),
        ("12", 122367809),
        ("12", 122445660),
        ("12", 124027929),
        ("12", 131933768),
        ("12", 131945789),
        ("13", 23141200),
        ("13", 30789743),
        ("14", 22899205),
        ("14", 23876742),
        ("14", 36205503),
        ("15", 23504535),
        ("15", 23513151),
        ("15", 75028597),
        ("15", 75131848),
        ("15", 83184713),
        ("15", 83925423),
        ("15", 83925828),
        ("15", 84085346),
        ("15", 88601453),
        ("16", 1328673),
        ("16", 2102840),
        ("16", 16199454),
        ("16", 23316941),
        ("16", 24073631),
        ("16", 45498322),
        ("16", 48911332),
        ("16", 67406951),
        ("17", 8157193),
        ("17", 17366356),
        ("17", 17661144),
        ("17", 30792312),
        ("17", 30792467),
        ("17", 35198045),
        ("17", 36981576),
        ("17", 55510846),
        ("17", 59254929),
        ("18", 2881847),
        ("18", 9172445),
        ("18", 9877857),
        ("18", 9878027),
        ("18", 61662156),
        ("19", 15713544),
        ("19", 41027164),
        ("19", 51499079),
        ("19", 51500424),
        ("20", 3076854),
        ("20", 20180356),
        ("20", 61796277),
        ("21", 39113508),
        ("22", 18338811),
        ("22", 18349495),
        ("22", 20328280),
        ("22", 20648538),
        ("22", 22912041),
        ("22", 29988205),
        ("22", 37425933),
        ("22", 43670607),
    )

    pth = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "data",
        "CEU.exon.2010_03.genotypes.trimmed.vcf.gz",
    )
    with open(pth, "rb") as vcf:
        pth = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "data",
            "CEU.exon.2010_03.genotypes.trimmed.vcf.gz.tbi",
        )
        with open(pth, "rb") as vcf_tbi:

            for _ in range(100):
                indexed = puretabix.TabixIndexedFile(vcf, vcf_tbi)

                # fetched = indexed.fetch("1", 1108138 - 10, 1108138 + 10)
                # fetched = indexed.fetch("1", 1108138)

                for chrom, pos in targets:
                    pos = int(pos)
                    fetched = indexed.fetch(chrom, pos)
                    print(fetched)